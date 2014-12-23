import pysam
from collections import Counter, namedtuple


def pos_counter_by_pos(bam_fpath, positions):
    alignmentfile = pysam.AlignmentFile(bam_fpath)
    ref_name_index = {}
    for pileup_col in alignmentfile.pileup():
        ref_id = pileup_col.reference_id
        ref_pos = pileup_col.reference_pos
        ref_name = ref_name_index.get(ref_id, alignmentfile.getrname(ref_id))
        position = (ref_name, ref_pos)
        if position in positions:
            pos_counter = Counter()
            print position
            for pileup_read in pileup_col.pileups:
                alignement = pileup_read.alignment
                blocks = alignement.get_blocks()

                if pileup_read.alignment.is_reverse:
                    base_pos = (pileup_read.alignment.query_length -
                                pileup_read.query_position - 1)
                else:
                    base_pos = pileup_read.query_position
#
#                 base_pos = pileup_read.query_position
                print base_pos, pileup_read.alignment.qname
                pos_counter[base_pos] += 1
            yield pos_counter


Block = namedtuple('Block', ['ref_start', 'ref_stop',
                             'read_start', 'read_stop'])


class ReadRefCoord(object):
    def __init__(self, alig_read, sam, hard_clip_as_soft=False):
        self._alig_read = alig_read
        self._sam = sam
        self.hard_clip_as_soft = hard_clip_as_soft
        self._blocks = None
        self._read_len = None
        self._last_position = None

    @property
    def read_len(self):
        if self._read_len is None:
            self.blocks
        return self._read_len

    @property
    def blocks(self):
        alig_read = self._alig_read
        if self._blocks is not None:
            return self._block
        ref_start = alig_read.reference_start

        cigar_tuples = alig_read.cigartuples
        is_fwd = not(alig_read.is_reverse)
        if not is_fwd:
            cigar_tuples = reversed(cigar_tuples)

        read_pos = 0  # if is reversed this is not true, but we fixed later
        read_len = 0
        ref_pos = ref_start
        blocks = []
        for operation, length in alig_read.cigartuples:
            blk_start = ref_pos, read_pos
            if operation == 5:
                if self.hard_clip_as_soft:
                    operation = 4   # if like a soft clip
                else:
                    continue
            if operation in (0, 7, 8):  # match
                read_len += length
                ref_pos += length
                if is_fwd:
                    read_pos += length
                else:
                    read_pos -= length
            elif operation == 1:  # insertion in read
                read_len += length
                if is_fwd:
                    read_pos += length
                else:
                    read_pos -= length
            elif operation == 2 or operation == 3:    # deletion in read
                ref_pos += length
            elif operation == 4:    # soft clip
                # This is not a block to be returned because it is not aligned
                read_len += length
                if is_fwd:
                    read_pos += length
                else:
                    read_pos -= length
            elif operation == 6:
                continue
            blk_end = ref_pos - 1, read_pos - 1
            if operation != 4:
                ref_start = blk_start[0]
                ref_end = blk_end[0]
                read_start = blk_start[1]
                read_end = blk_end[1]
                if operation in (2, 3):
                    read_start, read_end = None, None
                elif operation == 1:
                    ref_start, ref_end = None, None
                block = Block(ref_start, ref_end, read_start, read_end)
                blocks.append(block)

        if not is_fwd:
            rev_blocks = []
            for blk in blocks:
                if blk.read_stop is not None and blk.read_start is not None:
                    blk_read_start = blk.read_start + read_len - 1
                    blk_read_end = blk.read_stop + read_len + 1
                else:
                    blk_read_start = None
                    blk_read_end = None
                block = Block(blk.ref_start, blk.ref_stop, blk_read_start,
                              blk_read_end)
                rev_blocks.append(block)
            blocks = rev_blocks
        self._block = blocks
        self._read_len = read_len
        return blocks

    def get_read_pos(self, ref_pos):
        if (self._last_position is not None and
            self._last_position[0] == ref_pos):
            return self._last_position[1]
        alig_read = self._alig_read
        sam = self._sam
        chrom = sam.getrname(alig_read.reference_id)
        if chrom != ref_pos[0]:
            msg = 'The aligned read is not aligned to the given chrom: '
            msg += chrom
            raise ValueError(msg)
        ref_pos = ref_pos[1]

        blocks = self.blocks

        block_with_ref_pos = None
        for block in blocks:
            if block.ref_start <= ref_pos <= block.ref_stop:
                block_with_ref_pos = block
                break
        if block_with_ref_pos is None:
            read_pos = None
        else:
            read_start = block_with_ref_pos.read_start
            read_stop = block_with_ref_pos.read_stop

            if read_start is None or read_stop is None:
                read_pos = None
            else:
                ref_start = block_with_ref_pos.ref_start
                ref_stop = block_with_ref_pos.ref_stop

                if alig_read.is_reverse:
                    read_pos = read_start - (ref_pos - ref_start)
                else:
                    read_pos = read_stop - (ref_stop - ref_pos)

        self._last_position = ref_pos, read_pos
        return read_pos

    def get_read_pos_counting_from_end(self, ref_pos):
        if self._last_position is None:
            read_pos = self.get_read_pos(ref_pos)
        else:
            read_pos = self._last_position[1]
        read_len = self.read_len
        if read_pos is not None:
            return read_pos - read_len
